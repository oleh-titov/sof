// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright(c) 2021 Intel Corporation. All rights reserved.
//
// Author: Liam Girdwood <liam.r.girdwood@linux.intel.com>
// Author: Keyon Jie <yang.jie@linux.intel.com>

#include <sof/audio/buffer.h>
#include <sof/audio/component_ext.h>
#include <sof/audio/pipeline.h>
#include <sof/common.h>
#include <sof/drivers/idc.h>
#include <sof/ipc/topology.h>
#include <sof/ipc/common.h>
#include <sof/ipc/msg.h>
#include <sof/ipc/driver.h>
#include <sof/ipc/schedule.h>
#include <sof/lib/alloc.h>
#include <sof/lib/cache.h>
#include <sof/lib/cpu.h>
#include <sof/lib/mailbox.h>
#include <sof/list.h>
#include <sof/platform.h>
#include <sof/sof.h>
#include <sof/spinlock.h>
#include <ipc/dai.h>
#include <ipc/header.h>
#include <ipc/stream.h>
#include <ipc/topology.h>
#include <errno.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* create a new component in the pipeline */
struct comp_buffer *buffer_new(const struct sof_ipc_buffer *desc)
{
	struct comp_buffer *buffer;

	tr_info(&buffer_tr, "buffer new size 0x%x id %d.%d flags 0x%x",
		desc->size, desc->comp.pipeline_id, desc->comp.id, desc->flags);

	/* allocate buffer */
	buffer = buffer_alloc(desc->size, desc->caps, PLATFORM_DCACHE_ALIGN);
	if (buffer) {
		buffer->id = desc->comp.id;
		buffer->pipeline_id = desc->comp.pipeline_id;
		buffer->core = desc->comp.core;

		buffer->stream.underrun_permitted = desc->flags &
						    SOF_BUF_UNDERRUN_PERMITTED;
		buffer->stream.overrun_permitted = desc->flags &
						   SOF_BUF_OVERRUN_PERMITTED;

		memcpy_s(&buffer->tctx, sizeof(struct tr_ctx),
			 &buffer_tr, sizeof(struct tr_ctx));

		dcache_writeback_invalidate_region(buffer, sizeof(*buffer));
	}

	return buffer;
}

int32_t ipc_comp_pipe_id(const struct ipc_comp_dev *icd)
{
	switch (icd->type) {
	case COMP_TYPE_COMPONENT:
		return dev_comp_pipe_id(icd->cd);
	case COMP_TYPE_BUFFER:
		return icd->cb->pipeline_id;
	case COMP_TYPE_PIPELINE:
		return icd->pipeline->pipeline_id;
	default:
		tr_err(&ipc_tr, "Unknown ipc component type %u", icd->type);
		return -EINVAL;
	};

	return 0;
}

/* Function overwrites PCM parameters (frame_fmt, buffer_fmt, channels, rate)
 * with buffer parameters when specific flag is set.
 */
static void comp_update_params(uint32_t flag,
			       struct sof_ipc_stream_params *params,
			       struct comp_buffer *buffer)
{
	if (flag & BUFF_PARAMS_FRAME_FMT)
		params->frame_fmt = buffer->stream.frame_fmt;

	if (flag & BUFF_PARAMS_BUFFER_FMT)
		params->buffer_fmt = buffer->buffer_fmt;

	if (flag & BUFF_PARAMS_CHANNELS)
		params->channels = buffer->stream.channels;

	if (flag & BUFF_PARAMS_RATE)
		params->rate = buffer->stream.rate;
}

int comp_verify_params(struct comp_dev *dev, uint32_t flag,
		       struct sof_ipc_stream_params *params)
{
	struct list_item *buffer_list;
	struct list_item *source_list;
	struct list_item *sink_list;
	struct list_item *clist;
	struct list_item *curr;
	struct comp_buffer *sinkb;
	struct comp_buffer *buf;
	int dir = dev->direction;
	uint32_t flags = 0;

	if (!params) {
		comp_err(dev, "comp_verify_params(): !params");
		return -EINVAL;
	}

	source_list = comp_buffer_list(dev, PPL_DIR_UPSTREAM);
	sink_list = comp_buffer_list(dev, PPL_DIR_DOWNSTREAM);

	/* searching for endpoint component e.g. HOST, DETECT_TEST, which
	 * has only one sink or one source buffer.
	 */
	if (list_is_empty(source_list) != list_is_empty(sink_list)) {
		if (!list_is_empty(source_list))
			buf = list_first_item(&dev->bsource_list,
					      struct comp_buffer,
					      sink_list);
		else
			buf = list_first_item(&dev->bsink_list,
					      struct comp_buffer,
					      source_list);

		buffer_lock(buf, &flags);

		/* update specific pcm parameter with buffer parameter if
		 * specific flag is set.
		 */
		comp_update_params(flag, params, buf);

		/* overwrite buffer parameters with modified pcm
		 * parameters
		 */
		buffer_set_params(buf, params, BUFFER_UPDATE_FORCE);

		/* set component period frames */
		component_set_period_frames(dev, buf->stream.rate);

		buffer_unlock(buf, flags);
	} else {
		/* for other components we iterate over all downstream buffers
		 * (for playback) or upstream buffers (for capture).
		 */
		buffer_list = comp_buffer_list(dev, dir);
		clist = buffer_list->next;

		while (clist != buffer_list) {
			curr = clist;
			buf = buffer_from_list(curr, struct comp_buffer, dir);
			buffer_lock(buf, &flags);
			clist = clist->next;
			comp_update_params(flag, params, buf);
			buffer_set_params(buf, params, BUFFER_UPDATE_FORCE);
			buffer_unlock(buf, flags);
		}

		/* fetch sink buffer in order to calculate period frames */
		sinkb = list_first_item(&dev->bsink_list, struct comp_buffer,
					source_list);

		buffer_lock(sinkb, &flags);
		component_set_period_frames(dev, sinkb->stream.rate);
		buffer_unlock(sinkb, flags);
	}

	return 0;
}

int comp_buffer_connect(struct comp_dev *comp, uint32_t comp_core,
			struct comp_buffer *buffer, uint32_t dir)
{
	int ret;

	/* check if it's a connection between cores */
	if (buffer->core != comp_core) {
		dcache_invalidate_region(buffer, sizeof(*buffer));

		buffer->inter_core = true;

		if (!comp->is_shared) {
			comp = comp_make_shared(comp);
			if (!comp)
				return -ENOMEM;
		}
	}

	ret = pipeline_connect(comp, buffer, dir);
	dcache_writeback_invalidate_region(buffer, sizeof(*buffer));

	return ret;
}

int ipc_pipeline_complete(struct ipc *ipc, uint32_t comp_id)
{
	struct ipc_comp_dev *ipc_pipe;
	struct ipc_comp_dev *icd;
	struct pipeline *p;
	uint32_t pipeline_id;
	struct ipc_comp_dev *ipc_ppl_source;
	struct ipc_comp_dev *ipc_ppl_sink;
	int ret;

	/* check whether pipeline exists */
	ipc_pipe = ipc_get_comp_by_id(ipc, comp_id);
	if (!ipc_pipe) {
		tr_err(&ipc_tr, "ipc: ipc_pipeline_complete looking for pipe component id %d failed",
		       comp_id);
		return -EINVAL;
	}

	/* check core */
	if (!cpu_is_me(ipc_pipe->core))
		return ipc_process_on_core(ipc_pipe->core, false);

	p = ipc_pipe->pipeline;

	/* get pipeline source component */
	ipc_ppl_source = ipc_get_ppl_src_comp(ipc, p->pipeline_id);
	if (!ipc_ppl_source) {
		tr_err(&ipc_tr, "ipc: ipc_pipeline_complete looking for pipeline source failed");
		return -EINVAL;
	}

	/* get pipeline sink component */
	ipc_ppl_sink = ipc_get_ppl_sink_comp(ipc, p->pipeline_id);
	if (!ipc_ppl_sink) {
		tr_err(&ipc_tr, "ipc: ipc_pipeline_complete looking for pipeline sink failed");
		return -EINVAL;
	}

	/* find the scheduling component */
	icd = ipc_get_comp_by_id(ipc, p->sched_id);
	if (!icd) {
		tr_warn(&ipc_tr, "ipc_pipeline_complete(): no scheduling component specified, use comp %d",
			ipc_ppl_sink->id);

		icd = ipc_ppl_sink;
	}

	if (icd->type != COMP_TYPE_COMPONENT) {
		tr_err(&ipc_tr, "ipc_pipeline_complete(): icd->type (%d) != COMP_TYPE_COMPONENT for pipeline scheduling component icd->id %d",
		       icd->type, icd->id);
		return -EINVAL;
	}

	if (icd->core != ipc_pipe->core) {
		tr_err(&ipc_tr, "ipc_pipeline_complete(): icd->core (%d) != ipc_pipe->core (%d) for pipeline scheduling component icd->id %d",
		       icd->core, ipc_pipe->core, icd->id);
		return -EINVAL;
	}

	p->sched_comp = icd->cd;

	pipeline_id = ipc_pipe->pipeline->pipeline_id;

	tr_dbg(&ipc_tr, "ipc: pipe %d -> complete on comp %d", pipeline_id,
	       comp_id);

	ret = pipeline_complete(ipc_pipe->pipeline, ipc_ppl_source->cd,
				ipc_ppl_sink->cd);

	return ret;
}

int ipc_comp_free(struct ipc *ipc, uint32_t comp_id)
{
	struct ipc_comp_dev *icd;

	/* check whether component exists */
	icd = ipc_get_comp_by_id(ipc, comp_id);
	if (!icd)
		return -ENODEV;

	/* check core */
	if (!cpu_is_me(icd->core))
		return ipc_process_on_core(icd->core, false);

	/* check state */
	if (icd->cd->state != COMP_STATE_READY)
		return -EINVAL;

	/* free component and remove from list */
	comp_free(icd->cd);

	icd->cd = NULL;

	list_item_del(&icd->list);
	rfree(icd);

	return 0;
}
