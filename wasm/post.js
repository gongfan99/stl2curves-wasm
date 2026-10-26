// Keep lifetime management in JS: no exported C++ handles reach callers.
var stl2curvesDisposed = false;
var stl2curvesBusy = false;
var stl2curvesConvert = Module['convert'];
// Embind installs exports during initialization; wrap after runtime setup.
addOnPostRun(function () {
  stl2curvesConvert = Module['convert'];
  Module['convert'] = function (inputs, options) {
    if (stl2curvesDisposed) throw new Error('stl2curves module has been disposed');
    if (stl2curvesBusy) throw new Error('A conversion is already running in this module');
    stl2curvesBusy = true;
    try {
      return stl2curvesConvert(inputs, options === undefined ? {} : options);
    } finally {
      stl2curvesBusy = false;
    }
  };
});
Module['dispose'] = function () {
  if (stl2curvesDisposed) return;
  if (stl2curvesBusy) throw new Error('Terminate the application Worker to cancel an active conversion');
  stl2curvesDisposed = true;
  if (typeof PThread !== 'undefined') PThread.terminateAllThreads();
};
