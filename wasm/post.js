// Keep lifetime management in JS: no exported C++ handles reach callers.
var stl2curvesDisposed = false;
var stl2curvesBusy = false;
// In SDK 6.0.13, --post-js executes after the awaited runtime initialization.
// Pthread bootstrap instances do not expose the application Embind exports.
function stl2curvesGuard(name) {
  var operation = Module[name];
  if (typeof operation !== 'function') return;
  Module[name] = function () {
    if (stl2curvesDisposed) throw new Error('stl2curves module has been disposed');
    if (stl2curvesBusy) throw new Error('A model operation is already running in this module');
    stl2curvesBusy = true;
    try {
      var args = Array.prototype.slice.call(arguments);
      if ((name === 'convert' || name === 'exportGltf') && args[1] === undefined) args[1] = {};
      var result = operation.apply(null, args);
      if (result && result['__stl2curvesError']) {
        var error = result['__stl2curvesError'];
        throw error.name === 'TypeError' ? new TypeError(error.message) : new Error(error.message);
      }
      return result;
    } finally {
      stl2curvesBusy = false;
    }
  };
}
['convert', 'exportStep', 'exportGltf', 'getClosestPoints', 'getShapeProperties', 'releaseModel'].forEach(stl2curvesGuard);
Module['dispose'] = function () {
  if (stl2curvesDisposed) return;
  if (stl2curvesBusy) throw new Error('Terminate the application Worker to cancel an active model operation');
  stl2curvesDisposed = true;
  try {
    if (typeof Module['_clearModels'] === 'function') Module['_clearModels']();
  } finally {
    if (typeof PThread !== 'undefined') PThread.terminateAllThreads();
  }
};
