# To Do

* implement weighted sessions for Bayesian standardization
* implement drifts for Bayesian standardization

* Turn `D4xdata.R17_VPDB` (and perhaps also `D4xdata.R18_VPDB`) into a property, so that there is no risk of redefining `D4xdata.R17_VSMOW` without updating `D4xdata.R17_VPDB`; In this case, the setter methods for `D4xdata.R17_VPDB` and `D4xdata.R18_VPDB` should raise an exception.
* improve `test_virtual_data()` to populate with non-default parameters
* use a true CSV parser?
* add (many) plot customization options to CLI?
* implement individual weighing of replicates, e.g. according to the sqrt of the number of IRMS acquisitions
