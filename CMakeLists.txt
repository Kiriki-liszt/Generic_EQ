#- Check (target, target-au, project, elseif(SMTG_WIN) target_sources)

cmake_minimum_required(VERSION 3.14.0)

option(SMTG_ENABLE_VST3_PLUGIN_EXAMPLES "Enable VST 3 Plug-in Examples" OFF)
option(SMTG_ENABLE_VST3_HOSTING_EXAMPLES "Enable VST 3 Hosting Examples" OFF)

set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13 CACHE STRING "")

if(SMTG_MAC)
    set(vst3sdk_SOURCE_DIR "/Users/yoon/VST_SDK/vst3sdk")
elseif(SMTG_WIN)
    set(vst3sdk_SOURCE_DIR "C:/VST_SDK/vst3sdk")
elseif(GITHUB_ACTIONS) 
    set(vst3sdk_SOURCE_DIR ${CMAKE_SOURCE_DIR}/vst3sdk)
endif()

if(NOT vst3sdk_SOURCE_DIR)
    message(FATAL_ERROR "Path to VST3 SDK is empty!")
endif()

project(Generic_EQ
    # This is your plug-in version number. Change it here only.
    # Version number symbols usable in C++ can be found in
    # source/version.h and ${PROJECT_BINARY_DIR}/projectversion.h.
    VERSION 1.0.0.0
    DESCRIPTION "Generic EQ VST 3 Plug-in"
)

if(NOT SMTG_ENABLE_VSTGUI_SUPPORT)
    return()
endif()

set(target Generic_EQ)
set(target-au Generic_EQ-au)
set(bundleID Generic.EQ)

set(SMTG_VSTGUI_ROOT "${vst3sdk_SOURCE_DIR}")

add_subdirectory(${vst3sdk_SOURCE_DIR} ${PROJECT_BINARY_DIR}/vst3sdk)
smtg_enable_vst3_sdk()

smtg_add_vst3plugin(${target}
    source/version.h
    source/GNRC_EQ_cids.h
    source/GNRC_EQ_shared.h        # added
    source/GNRC_EQ_processor.h
    source/GNRC_EQ_processor.cpp
    source/GNRC_EQ_controller.h
    source/GNRC_EQ_controller.cpp
    source/GNRC_EQ_entry.cpp
    resource/GNRC_EQ_editor.uidesc
)

smtg_target_add_plugin_resources(${target}
    RESOURCES
        resource/GNRC_EQ_editor.uidesc
        # resource/background.png
)

smtg_target_add_plugin_snapshots (${target}
    RESOURCES
        resource/BE9ACDC3EA6D5D86A6A76A31FAEED6B4_snapshot.png
        resource/BE9ACDC3EA6D5D86A6A76A31FAEED6B4_snapshot_2.0x.png
)

target_compile_features(${target}
    PUBLIC
        cxx_std_20 
)

target_link_libraries(${target}
    PRIVATE
        sdk
        vstgui_support
)

smtg_target_configure_version_file(${target})

if(SMTG_MAC)
    smtg_target_set_bundle(${target}
        BUNDLE_IDENTIFIER "io.github.yg331.${bundleID}"
        COMPANY_NAME "yg331"
    )
    smtg_target_set_debug_executable(${target}
        "/Applications/VST3PluginTestHost.app"
        "--pluginfolder;$(BUILT_PRODUCTS_DIR)"
    )
elseif(SMTG_WIN)
    target_sources(${target} 
        PRIVATE 
            resource/win32resource.rc
            resource/GNRC_EQ.rc
    )
    if(MSVC)
        set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${target})

        smtg_target_set_debug_executable(${target}
            "$(ProgramW6432)/Steinberg/VST3PluginTestHost/VST3PluginTestHost.exe"
            "--pluginfolder \"$(OutDir)/\""
        )
    endif()
endif(SMTG_MAC)

# Add an AUv2 target
if (SMTG_MAC AND XCODE AND SMTG_ENABLE_AUV2_BUILDS)
    list(APPEND CMAKE_MODULE_PATH "${vst3sdk_SOURCE_DIR}/cmake/modules")
    include(SMTG_AddVST3AuV2) 
	smtg_target_add_auv2(${target-au}
		BUNDLE_NAME ${target}
		BUNDLE_IDENTIFIER io.github.yg331.${bundleID}.audiounit
		INFO_PLIST_TEMPLATE resource/au-info.plist
		VST3_PLUGIN_TARGET ${target})
endif(SMTG_MAC AND XCODE AND SMTG_ENABLE_AUV2_BUILDS)

