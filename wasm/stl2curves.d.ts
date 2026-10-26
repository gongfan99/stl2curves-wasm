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
export interface ConversionResult {
  files: OutputFile[];
  errors: { input: string; object: string; code: string; message: string; diagnostics?: Diagnostics }[];
}
export interface Stl2curvesModule {
  convert(inputs: InputFile[], options?: ConversionOptions): ConversionResult;
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
