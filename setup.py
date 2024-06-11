import setuptools
import os


setuptools.setup(

	name="pv3rs_wrapper",

	version="0.0.1",
	packages=["pv3rs_wrapper"],
	license="Pv3Rs",
	long_description="Pv3Rs wrapper command line tool",
	scripts= ["scripts/%s" % x for x in os.listdir("scripts")],
)
