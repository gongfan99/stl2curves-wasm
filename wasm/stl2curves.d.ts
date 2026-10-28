export interface InputFile { name: string; data: Uint8Array }
export interface ConversionOptions {
  tolerance?: number;
  accuracy?: number;
  timeLimit?: number;
  /** -1 automatic, 0 disabled, positive values specify millimeters. */
  simplify?: number;
  fuse?: boolean;
  curves?: boolean;
  blends?: boolean;
  repair?: boolean;
  trueSize?: boolean;
  details?: boolean;
  merge?: string;
}
export interface Diagnostics {
  inputTriangles: number; triangles: number; faces: number; solids: number;
  openEdges: number; nonmanifoldEdges: number; parallelTasks: number;
  meshVolume: number; volume: number; scale: number;
  designUnit: '' | 'mm' | 'inch'; designLabel: string;
  valid: boolean; fallback: boolean; timeLimited: boolean;
  repairs: string[]; warnings: string[]; fallbacks: string[];
  surfaceCounts: Record<string, number>;
  features: { kind: string; triangles: number; deviation: number; radius: number; pitch: number }[];
}
export interface OutputFile { name: string; data: Uint8Array<ArrayBuffer>; diagnostics: Diagnostics }
export interface GraphId { kind: number; index: number }
export interface UsagePathStep { node: GraphId; ref?: GraphId; stepIndex: number }
export interface ShapeRef { modelId: number; brepgraph_id?: GraphId; usagePath?: UsagePathStep[] }
export interface ModelInfo { modelId: number; name: string; root: ShapeRef; diagnostics: Diagnostics }
export interface ConversionError { input: string; object: string; code: string; message: string; diagnostics?: Diagnostics }
export interface GltfExportOptions {
  /** Maximum linear deflection in model millimeters; default 0.1. */
  linearDeflection?: number;
  /** Angular deflection in radians; default 0.5. */
  angularDeflection?: number;
  faces?: boolean; edges?: boolean; vertices?: boolean;
}
export interface GltfFile { name: string; data: Uint8Array<ArrayBuffer>; warnings: string[] }
export interface ConversionResult {
  models: ModelInfo[];
  errors: ConversionError[];
}
export interface Stl2curvesModule {
  convert(inputs: InputFile[], options?: ConversionOptions): ConversionResult;
  exportStep(modelId: number): OutputFile;
  exportGltf(modelId: number, options?: GltfExportOptions): GltfFile;
  getClosestPoints(a: ShapeRef, b: ShapeRef): { p1: [number, number, number]; p2: [number, number, number]; distance: number };
  getShapeProperties(shape: ShapeRef): { properties: { name: string; value: number[] }[] };
  releaseModel(modelId: number): boolean;
  buildInfo(): { version: string; occtVersion: string; threaded: boolean; workers: number };
  dispose(): void;
}
export interface ModuleConfig {
  locateFile?: (path: string, prefix: string) => string;
  mainScriptUrlOrBlob?: string | Blob;
  wasmBinary?: ArrayBuffer | Uint8Array;
  print?: (message: string) => void;
  printErr?: (message: string) => void;
}
export default function createStl2curves(config?: ModuleConfig): Promise<Stl2curvesModule>;
