workspace "main"
    architecture "x64"
    toolset "gcc"

    buildoptions "--std=c++23"

    configurations {
        "debug",
        "release",
        "dist"
    }

    filter "configurations:debug"
        defines {"_DEBUG", "_PROFILING"}
        symbols "On"
        filter {"system:linux", "action:gmake"}
            buildoptions {"-g"} -- For core files
    filter "configurations:release"
        defines {"_RELEASE", "NDEBUG", "_PROFILING"}
        optimize "On"
    filter "configurations:dist"
        defines {"_DIST", "NDEBUG", "_PROFILING"}
        optimize "Speed"

    filter "system:macosx"
        defines {"_MACOSX"}

    filter "system:linux"
        defines {"_LINUX"}

    filter "system:windows"
        defines {"_WINDOWS"}

outputdir = "%{cfg.buildcfg}/%{cfg.system}/%{cfg.architectrure}"
project "VMesh"
    kind "StaticLib"
    language "C++"
    targetname "VMesh"
    targetdir ("bin/" .. outputdir)
    objdir ("bin-int/" .. outputdir)

    files {
        "src/**.hpp",
        "src/**.cpp",
        "include/**.hpp",
        "include/**.cpp",
        "deps/src/**.cpp",
        "deps/src/**.c",
        "deps/include/**.h",
        "deps/include/**.hpp"
    }

    includedirs {
        "include",
        "deps/include",
        "/usr/include"
    }

    links {
        "assimp"
    }

