import setuptools
import os


setuptools.setup(

	name="DEploid_Pv_pipeline",

	version="0.0.1",
	packages=["DEploid_Pv"],
	license="DEploid",
	long_description="DEploid wrapper command line tool",
	scripts= ["scripts/%s" % x for x in os.listdir("scripts")],
)
