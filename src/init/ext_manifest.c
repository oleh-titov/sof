// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright(c) 2020 Intel Corporation. All rights reserved.
//
// Author: Karol Trzcinski <karolx.trzcinski@linux.intel.com>
//

#include <sof/bit.h>
#include <sof/common.h>
#include <sof/compiler_info.h>
#include <sof/debug/debug.h>
#include <kernel/abi.h>
#include <kernel/ext_manifest.h>
#include <config.h>
#include <version.h>

const struct ext_man_fw_version ext_man_fw_ver
	__aligned(EXT_MAN_ALIGN) __section(".fw_metadata") = {
	.hdr.type = EXT_MAN_ELEM_FW_VERSION,
	.hdr.elem_size = ALIGN_UP(sizeof(struct ext_man_fw_version),
				  EXT_MAN_ALIGN),
	.version = {
		.hdr.size = sizeof(struct sof_ipc_fw_version),
		.micro = SOF_MICRO,
		.minor = SOF_MINOR,
		.major = SOF_MAJOR,
#if CONFIG_DEBUG
		/* only added in debug for reproducibility in releases */
		.build = SOF_BUILD,
		.date = __DATE__,
		.time = __TIME__,
#endif
		.tag = SOF_TAG,
		.abi_version = SOF_ABI_VERSION,
	},
	.flags = DEBUG_SET_FW_READY_FLAGS,
};

const struct ext_man_cc_version ext_man_cc_ver
	__aligned(EXT_MAN_ALIGN) __section(".fw_metadata") = {
	.hdr.type = EXT_MAN_ELEM_CC_VERSION,
	.hdr.elem_size = ALIGN_UP(sizeof(struct ext_man_cc_version),
				  EXT_MAN_ALIGN),
	.cc_version = {
		.ext_hdr.hdr.size = sizeof(struct sof_ipc_cc_version),
		.ext_hdr.hdr.cmd = SOF_IPC_FW_READY,
		.ext_hdr.type = SOF_IPC_EXT_CC_INFO,
		.micro = CC_MICRO,
		.minor = CC_MINOR,
		.major = CC_MAJOR,
		.name = CC_NAME "\0",	///< eg. "XCC", "\0" is needed when
					///< sizeof(CC_NAME)-1 == sizeof(.name)
		.optim = CC_OPTIMIZE_FLAGS "\0", ///< eg. "O2"
		.desc = CC_DESC "\0", ///< eg. " RG-2017.8-linux"
	},
};
