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
      "haystackfm-0.5.0" = "sha256-yKWwEjemqm+/RCufg3aRVAgfpfg3vs+biDwdcKgwri4=";
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
