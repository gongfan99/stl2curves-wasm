param(
    [ValidateSet('single_thread', 'threaded', 'all')]
    [string]$Variant = 'single_thread',
    [string]$Workspace = 'E:\wasm',
    [ValidateRange(1, 64)]
    [int]$Jobs = 2
)
$ErrorActionPreference = 'Stop'
$configPath = Join-Path $Workspace 'local-build.json'
if (-not (Test-Path -LiteralPath $configPath)) {
    throw "Missing local build configuration: $configPath"
}
$config = Get-Content -LiteralPath $configPath -Raw | ConvertFrom-Json
$repository = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot '..')).Path
$linuxRepository = (& wsl.exe -d $config.distribution --exec wslpath -a $repository).Trim()
if ($LASTEXITCODE -ne 0) { throw 'Cannot resolve the source repository in WSL' }
$linuxOutput = (& wsl.exe -d $config.distribution --exec wslpath -w "$($config.linuxWorkspace)/dist").Trim()
if ($LASTEXITCODE -ne 0) { throw 'Cannot resolve the WASM output directory in WSL' }
& wsl.exe -d $config.distribution --cd $linuxRepository --exec bash scripts/build_wasm_local.sh $Variant $config.linuxWorkspace $Jobs
if ($LASTEXITCODE -ne 0) { throw "WASM build failed with exit code $LASTEXITCODE" }
$destination = Join-Path $Workspace 'dist'
if (-not [string]::Equals($linuxOutput, $destination, [StringComparison]::OrdinalIgnoreCase)) {
    $source = '\\wsl.localhost\' + $config.distribution + ($config.linuxWorkspace -replace '/', '\') + '\dist'
    New-Item -ItemType Directory -Force -Path $destination | Out-Null
    foreach ($entry in Get-ChildItem -LiteralPath $source) {
        Copy-Item -LiteralPath $entry.FullName -Destination $destination -Recurse -Force
    }
}
Write-Host "WASM artifacts: $destination"
