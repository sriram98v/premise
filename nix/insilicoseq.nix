# InSilicoSeq — the read simulator behind every synthetic sample in the benchmark.
#
{ lib
, python3Packages
, fetchPypi
}:

python3Packages.buildPythonApplication rec {
  pname = "insilicoseq";
  version = "2.0.1";
  pyproject = true;

  src = fetchPypi {
    pname = "InSilicoSeq";
    inherit version;
    hash = "sha256-59ZJP1+tTeBZ01WncEHZsHfJ7cE+wDVSrqpPUntvTpE=";
  };

  build-system = [ python3Packages.setuptools ];

  dependencies = with python3Packages; [
    numpy
    scipy
    biopython
    pysam
    requests
  ];

  doCheck = false;
  pythonImportsCheck = [ "iss" "iss.generator" ];

  postInstall = ''
    profile=$(find $out -name 'miSeq_0.npz' -print -quit)
    if [ -z "$profile" ]; then
      echo "insilicoseq: miSeq_0.npz missing from the install — the benchmark's" >&2
      echo "  --model miseq generators would fail at runtime." >&2
      exit 1
    fi
    echo "insilicoseq: error model present at $profile"
  '';

  doInstallCheck = true;
  postInstallCheck = ''
    $out/bin/iss --version
  '';

  meta = with lib; {
    description = "Sequencing read simulator producing realistic Illumina reads";
    homepage = "https://github.com/HadrienG2/InSilicoSeq";
    license = licenses.mit;
    mainProgram = "iss";
    platforms = platforms.unix;
  };
}
