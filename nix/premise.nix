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
  };

  buildType = "release";

  doCheck = false;

  meta = with lib; {
    description = "Paired-end metagenomic read classifier (seed -> align -> EM)";
    homepage = "https://github.com/sriram98v/premise";
    license = licenses.mit;
    mainProgram = "premise";
    platforms = platforms.unix;
  };
}
