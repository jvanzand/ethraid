# I don't know why, but I have to do 
#    from some_module import *
# in order for the kern_profiler to work. In particular, the @profile decorators throw an error when I run setup.py UNLESS I have the import statement above. I need to get to the bottom of this at some point, but for now it seems harmless.