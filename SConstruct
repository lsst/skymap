# -*- python -*-
#
# Setup our environment
#
import os
import lsst.SConsUtils as scons

try:
    scons.ConfigureDependentProducts
except AttributeError:
    import lsst.afw.SconsUtils
    scons.ConfigureDependentProducts = lsst.afw.SconsUtils.ConfigureDependentProducts

env = scons.makeEnv("skymap",
                    r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/skymap/trunk/SConstruct $",
                    scons.ConfigureDependentProducts("skymap"))

env.libs["skymap"] += env.getlibs("boost wcslib gsl utils daf_base daf_data") + \
    env.getlibs("daf_persistence pex_exceptions pex_logging pex_policy security afw")

env['IgnoreFiles'] = r"(~$|\.pyc$|^\.svn$|\.o$)"

# Locate SConscript files
for d in Split("doc examples lib python/lsst/skymap tests"):
    sconscriptPath = os.path.join(d, "SConscript")
    if os.path.exists(sconscriptPath):
        SConscript(os.path.join(d, "SConscript"))

# Install all directories
specialDirList = ("doc", "ups") # need special installation
stdDirList = [dirName for dirName in os.listdir(".") if os.path.isdir(dirName)
    and dirName not in specialDirList and not dirName.startswith(".")]
Alias("install", [
    env.InstallAs(os.path.join(env['prefix'], "doc", "doxygen"), os.path.join("doc", "htmlDir")),
    env.InstallEups(env['prefix'] + "/ups"),
] + [env.Install(env['prefix'], stdDir) for stdDir in stdDirList])

# What files should be cleaned
scons.CleanTree(r"*~ core *.so *.os *.o *.pyc")

files = scons.filesToTag()
if files:
    env.Command("TAGS", files, "etags -o $TARGET $SOURCES")

env.Declare()
env.Help("""
Store and retrieve all-sky image data
""")
