# fit-bpass-models
Fit BPASS sed models to observed photometric data.

For the time being, keep the binary fraction in the SSP fixed.

The free parameters are then age, metallicity, and escape fraction.

Currently, escape fraction is implemented as affecting F336W alone.

Need to include lines as well.

Try using the escape fraction to set the line strength and the F336W strength. Try using the input SED lyman continuum photon production rate to power the lines (x 1-fesc).