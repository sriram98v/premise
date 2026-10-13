{ lib
, rustPlatform
, src
, version
}:

rustPlatform.buildRustPackage {
  pname = "premise";
  inherit version src;

  cargoLock = {
    lockFile = "${src}/Cargo.lock";
    outputHashes = {
      "haystackfm-0.6.8" = "sha256-1knSJSCm1SK0LTkH3JiCMyZ5Jo25rjBGxl7YR26ryYQ=";
    };
  };

  buildType = "release";

  doCheck = false;

  meta = with lib; {
    description = "PREMISE: A probabilistic framework for source assignment of Illumina reads";
    homepage = "https://github.com/sriram98v/premise";
    license = licenses.mit;
    mainProgram = "premise";
    platforms = platforms.unix;
  };
}
