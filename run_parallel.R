# .libPaths("/data/nordicmathcovid/rpackages/")
# Uncomment to run targets sequentially on your local machine.
targets::tar_make_future(workers = 8, reporter = "timestamp")

q(save = "no")