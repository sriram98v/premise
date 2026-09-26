# Derivations for the compared methods.

{ lib, stdenv, fetchFromGitHub, fetchurl, fetchgit, cmake, zlib, bzip2, yaml-cpp, python3Packages, rustPlatform
, autoPatchelfHook, ... }:

let
  xxHashSrc = fetchFromGitHub {
    owner = "Cyan4973";
    repo = "xxHash";
    rev = "f9155bd4c57e2270a4ffbb176485e5d713de1c9b";
    hash = "sha256-YizwYMS2BLQlakCGpXX55F9bXinGO5rCepqSDc6EhpY=";
  };

  fixSeqan3Odr = seqan3RelPath: ''
    if [ -d "${seqan3RelPath}/include/seqan3" ]; then
      find "${seqan3RelPath}/include/seqan3" -name '*.hpp' -print0 \
        | xargs -0 sed -i -E 's/(^|[^e] )constexpr bool add_enum_bitwise_operators</\1inline constexpr bool add_enum_bitwise_operators</'
      echo "patched seqan3 ODR specialisations under ${seqan3RelPath}"
    fi
  '';
in
rec {
  raptor = stdenv.mkDerivation rec {
    pname = "raptor";
    version = "3.0.1";
    src = fetchurl {
      url = "https://github.com/seqan/raptor/releases/download/raptor-v${version}/raptor-${version}-Source.tar.xz";
      hash = "sha256-vkTnsmNeEQyUDPQP8TKVYgBVGi4ldl2fBX6ET3o8zrs=";
    };

    nativeBuildInputs = [ cmake ];
    buildInputs = [ zlib yaml-cpp ];
    enableParallelBuilding = true;

    postPatch = fixSeqan3Odr "lib/seqan3" + ''
      substituteInPlace lib/chopper/CMakeLists.txt \
        --replace-fail 'GIT_REPOSITORY https://github.com/Cyan4973/xxHash.git' \
                       'SOURCE_DIR "${xxHashSrc}"' \
        --replace-fail 'GIT_TAG "f9155bd4c57e2270a4ffbb176485e5d713de1c9b"' ""
    '';

    cmakeFlags = [
      "-DCMAKE_POLICY_VERSION_MINIMUM=3.5"
    ];

    env.CXXFLAGS = "-DSEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC -include cstdint";
    env.CFLAGS = "-include stdint.h";

    env.CMAKE_POLICY_VERSION_MINIMUM = "3.5";

    buildFlags = [ "raptor" ];

    installPhase = ''
      runHook preInstall
      mkdir -p $out/bin
      install -m755 bin/raptor $out/bin/raptor
      runHook postInstall
    '';

    meta = {
      description = "ganon's HIBF backend";
      homepage = "https://github.com/seqan/raptor";
      license = lib.licenses.bsd3;
      platforms = lib.platforms.unix;
      mainProgram = "raptor";
      };
  };

  ganon-cpp = stdenv.mkDerivation rec {
    pname = "ganon-cpp";
    version = "2.4.2";
    src = fetchgit {
      url = "https://github.com/pirovc/ganon.git";
      rev = "v${version}";
      fetchSubmodules = true;
      hash = "sha256-t37TvAHhDZJM40aEMzldZfYmHlBbbniQn7i0T01jTW8=";
    };
    nativeBuildInputs = [ cmake ];

    buildInputs = [ zlib bzip2 ];
    enableParallelBuilding = true;

    postPatch = fixSeqan3Odr "libs/seqan3" + ''
      # Write the version explicitly from the pinned tag since default falls to 0.0.0.
      echo "${version}" > VERSION.txt
    '';

    env.CXXFLAGS = "-DSEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC";
    buildFlags = [ "ganon-build" "ganon-classify" ];

    installPhase = ''
      runHook preInstall
      mkdir -p $out/bin
      install -m755 ganon-build ganon-classify $out/bin/
      runHook postInstall
    '';
  };

  pylowestcommonancestor = python3Packages.buildPythonPackage rec {
    pname = "pylowestcommonancestor";
    version = "1.0.0";
    pyproject = true;
    src = python3Packages.fetchPypi {
      inherit pname version;
      hash = "sha256-OExwB6AEcqRC0I1YXWPRp8LG2CJSBI5iPyc+9TuTD+Y=";
    };
    build-system = [ python3Packages.setuptools ];
    doCheck = false;
  };

  multitax = python3Packages.buildPythonPackage rec {
    pname = "multitax";
    version = "1.6.0";
    pyproject = true;
    src = python3Packages.fetchPypi {
      inherit pname version;
      hash = "sha256-Y0Evc2Y+le+mxOF06eMBr+zCHT85PaJgZX120f3oHqQ=";
    };
    build-system = [ python3Packages.setuptools python3Packages.setuptools-scm ];
    dependencies = [ pylowestcommonancestor ];
    doCheck = false;
    pythonImportsCheck = [ "multitax" ];
  };

  ganon = python3Packages.buildPythonApplication rec {
    pname = "ganon";
    version = "2.4.2";
    pyproject = true;
    src = ganon-cpp.src;
    build-system = [ python3Packages.setuptools python3Packages.setuptools-scm ];
    dependencies = [ python3Packages.pandas multitax ];

    SETUPTOOLS_SCM_PRETEND_VERSION_FOR_GANON = version;
    doCheck = false;
    makeWrapperArgs = [ "--prefix" "PATH" ":" "${ganon-cpp}/bin" "--prefix" "PATH" ":" "${raptor}/bin" ];
    postInstallCheck = ''
      $out/bin/ganon --version
    '';
  };

  karp = stdenv.mkDerivation {
    pname = "karp";
    version = "0-unstable-2017-09-16";
    src = fetchFromGitHub {
      owner = "mreppell";
      repo = "Karp";
      rev = "88c5b14";
      hash = "sha256-ibvRHgVblJemHQJ0QCyGWCRX9Wjcb7EIlW7ekBz+Ilk=";
    };
    nativeBuildInputs = [ cmake ];
    buildInputs = [ zlib ];
    enableParallelBuilding = true;

    # from claude
    # Two missing includes, without which karp does not compile at all against a current libstdc++:
    #   fastaIndex.h     needs <cstdint> — otherwise uint64_t is undeclared, findex_entry's
    #                                     `linelength` member never comes into existence, and every
    #                                     use site in ProcessReads.cpp fails with
    #                                     "'...findex_entry' has no member named 'linelength'"
    #   MinCollector.cpp needs <limits>
    # These are not upstream and were recorded nowhere: they existed only as uncommitted edits in
    # tools/karp/, a directory that is not under version control. The karp binary behind the published
    # tables was built from them, which means those results were not reproducible from any pinned
    # source until this patch captured the edits.
    patches = [ ./patches/karp-missing-includes.patch ];

    # karp declares find_package(HDF5 REQUIRED) and links ${HDF5_LIBRARIES} into both targets, but no
    # karp source includes an HDF5 header — inherited from kallisto, which karp derives from. Rather
    # than supply a library it never calls, the requirement is satisfied with empty lists. (Supplying
    # nixpkgs' hdf5 does not work: module-mode FindHDF5 picks an h5cc wrapper that cannot see across
    # nix's split outputs, and config mode hard-errors on HDF5_TOOLS_DIR = <out>/bin, which does not
    # exist because the tools live in a separate output.)
    postPatch = ''
      substituteInPlace src/CMakeLists.txt \
        --replace 'find_package( HDF5 REQUIRED )' \
                  'set( HDF5_FOUND TRUE )
    set( HDF5_LIBRARIES "" )
    set( HDF5_INCLUDE_DIR "" )'
    '';

    cmakeFlags = [ "-DCMAKE_POLICY_VERSION_MINIMUM=3.5" ];

    env.NIX_CFLAGS_LINK = "-lz";

    installPhase = ''
      runHook preInstall
      mkdir -p $out/bin
      install -m755 src/karp $out/bin/karp
      runHook postInstall
    '';
  };

  centrifuger = stdenv.mkDerivation {
    pname = "centrifuger";
    version = "1.1.3-r347";
    src = fetchFromGitHub {
      owner = "mourisl";
      repo = "centrifuger";
      rev = "ff9b4f41f5ffe4163fe682833337376ad8fd81ee";
      hash = "sha256-ND5DoQhDOtXT2mq97AS0k2yie19CLaX1b4V6D4IewzM=";
    };
    buildInputs = [ zlib ];
    enableParallelBuilding = true;
    installPhase = ''
      runHook preInstall
      mkdir -p $out/bin
      install -m755 centrifuger centrifuger-build centrifuger-inspect centrifuger-download $out/bin/ 2>/dev/null || \
        install -m755 centrifuger centrifuger-build $out/bin/
      runHook postInstall
    '';
  };

  mora = rustPlatform.buildRustPackage {
    pname = "mora";
    version = "1.0.0";
    src = fetchFromGitHub {
      owner = "AfZheng126";
      repo = "MORA";
      rev = "c06fa3d9d04d0705edfabc2948cf93a2ff742f51";
      hash = "sha256-6T8O2sDJP3UzmCjncDnlpvscNbeGJuZjzg5jnPZx+4U=";
    };

    nativeBuildInputs = [ cmake ];
    cargoHash = "sha256-tT1vf7aC3ZvvEI+52VCbbkVx3qKvgctdpOwHEBGpO+U=";
    buildInputs = [ zlib ];

    LIBZ_SYS_STATIC = "0";
    doCheck = false;
  };

  sylph = rustPlatform.buildRustPackage rec {
    pname = "sylph";
    version = "0.9.0";
    src = fetchFromGitHub {
      owner = "bluenote-1577";
      repo = "sylph";
      rev = "v${version}";
      hash = "sha256-HeHl5Oe8qcRRqnyWKhWJxZ7LtLPQHGkpDv890YXD8Jo=";
    };

    nativeBuildInputs = [ cmake ];
    cargoHash = "sha256-xY2Cj43np9J0ayJ90RNsA7aPah87MOhkuiJRZBe2PsE=";
    doCheck = false;
  };

  kmcp = stdenv.mkDerivation rec {
    pname = "kmcp";
    version = "0.9.5";
    src = fetchurl {
      url = "https://github.com/shenwei356/kmcp/releases/download/v${version}/kmcp_linux_amd64.tar.gz";
      hash = "sha256-0caWmeArXjiTIG9+oVbMNb5AZYm8vs44l9Ad03netFQ=";
    };
    sourceRoot = ".";
    nativeBuildInputs = [ autoPatchelfHook ];
    installPhase = ''
      runHook preInstall
      mkdir -p $out/bin
      install -m755 kmcp $out/bin/kmcp
      runHook postInstall
    '';
  };

  salmon = stdenv.mkDerivation rec {
    pname = "salmon";
    version = "2.5.1";
    src = fetchurl {
      url = "https://github.com/COMBINE-lab/salmon/releases/download/v${version}/"
            + "salmon-cli-x86_64-unknown-linux-gnu.tar.xz";
      hash = "sha256-atKgGyAiCS+I9MlWAdGXUeWU3TqnvnfufBHBq/SELL4=";
    };
    nativeBuildInputs = [ autoPatchelfHook ];
    buildInputs = [ stdenv.cc.cc.lib ];
    installPhase = ''
      runHook preInstall
      mkdir -p $out/bin
      install -m755 salmon $out/bin/salmon
      runHook postInstall
    '';
  };
}
