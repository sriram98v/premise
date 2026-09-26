{
  #
  #     nix develop                    # dev shell: build and hack on premise
  #     nix develop .#benchmark        # study shell: every tool `python3 -m premise_bench` (comparative-analysis/) needs
  #
  description = "premise — paired-end metagenomic read classifier, and its benchmark toolchain";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      lib = nixpkgs.lib;

      pkgs = import nixpkgs {
        inherit system;
        config.allowUnfreePredicate = p:
          builtins.elem (lib.getName p) [ "sratoolkit" "ncbi-vdb" "corefonts" ];
      };
      py = pkgs.python3Packages;

      methods = pkgs.callPackage ./nix/methods.nix { stdenv = pkgs.gcc13Stdenv; };
      insilicoseq = pkgs.callPackage ./nix/insilicoseq.nix { };

      premise-pkg = pkgs.callPackage ./nix/premise.nix {
        src = self;
        version = "0.4.0";
      };

      xopen = py.buildPythonPackage rec {
        pname = "xopen";
        version = "2.1.0";
        pyproject = true;
        src = pkgs.fetchPypi {
          inherit pname version;
          hash = "sha256-BlhIIViPgTVjhjtj1eAlaKQyhmG3UH6jjEMvJMaepLc=";
        };
        build-system = [ py.setuptools py.setuptools-scm ];
        dependencies = [ py.isal py.zlib-ng ];
        doCheck = false;
        pythonImportsCheck = [ "xopen" ];
      };

      dnaio = py.buildPythonPackage rec {
        pname = "dnaio";
        version = "1.2.4";
        pyproject = true;
        src = pkgs.fetchPypi {
          inherit pname version;
          hash = "sha256-p1cDEfKeizweo5pg9Xt7r42tjyUIWVxY1CeMVXFGMWY=";
        };
        build-system = [ py.setuptools py.setuptools-scm py.cython ];
        dependencies = [ xopen ];
        doCheck = false;
        pythonImportsCheck = [ "dnaio" ];
      };

      fastapy = py.buildPythonPackage rec {
        pname = "fastapy";
        version = "1.0.5";
        pyproject = true;
        src = pkgs.fetchPypi {
          inherit pname version;
          hash = "sha256-Ve92sN7ZBxQgAj1kylL4V7iR24UVev4+i+ew9+oCfxs=";
        };
        build-system = [ py.setuptools py.wheel ];
        doCheck = false;
        pythonImportsCheck = [ "fastapy" ];
      };

      cutadapt = py.buildPythonApplication rec {
        pname = "cutadapt";
        version = "5.2";
        pyproject = true;
        src = pkgs.fetchPypi {
          inherit pname version;
          hash = "sha256-I5Te6tQuyuX+D982njXz4q/tdw4UBZWCJyd5wugpXTw=";
        };
        build-system = [ py.setuptools py.setuptools-scm py.cython ];
        dependencies = [ dnaio xopen ];
        doCheck = false;

        postInstallCheck = ''
          $out/bin/cutadapt --version
        '';
      };

      bwa-mem2 = pkgs.gcc13Stdenv.mkDerivation rec {
        pname = "bwa-mem2";
        version = "2.3";

        src = pkgs.fetchFromGitHub {
          owner = "bwa-mem2";
          repo = "bwa-mem2";
          rev = "v${version}";
          fetchSubmodules = true;
          hash = "sha256-Oq3VBeFGrTa6hUUArGJcX+WECM8jiPtzohh/2L5VPu0=";
        };

        buildInputs = [ pkgs.zlib ];
        enableParallelBuilding = true;
        installPhase = ''
          runHook preInstall
          mkdir -p $out/bin
          for b in bwa-mem2 bwa-mem2.sse41 bwa-mem2.sse42 bwa-mem2.avx bwa-mem2.avx2 bwa-mem2.avx512bw; do
            [ -f "$b" ] && install -m755 "$b" $out/bin/
          done
          test -x $out/bin/bwa-mem2
          runHook postInstall
        '';
      };

      benchPython = pkgs.python3.withPackages (ps: [
        ps.matplotlib
        ps.pandas
        ps.pysam
        fastapy
        (py.toPythonModule insilicoseq)
      ]);

      benchFontsConf = pkgs.makeFontsConf { fontDirectories = [ pkgs.corefonts ]; };

      benchmarkBins = [
        premise-pkg

        methods.karp
        methods.centrifuger
        methods.mora
        methods.sylph
        methods.kmcp
        methods.ganon
        methods.ganon-cpp
        methods.raptor
        methods.salmon

        bwa-mem2
        cutadapt
        pkgs.sratoolkit
        pkgs.samtools
        pkgs.seqkit

        pkgs.corefonts
      ];
    in {
      packages.${system} = {
        inherit cutadapt dnaio xopen fastapy bwa-mem2 insilicoseq;
        inherit (methods) raptor ganon ganon-cpp multitax pylowestcommonancestor karp centrifuger mora sylph kmcp salmon;
        premise = premise-pkg;
        default = premise-pkg;

        benchmarkTools = pkgs.buildEnv {
          name = "premise-benchmark-tools";
          paths = benchmarkBins ++ [ benchPython ];
        };
      };

      devShells.${system} = {

        default = pkgs.mkShell {
          packages = [
            pkgs.cargo
            pkgs.rustc
            pkgs.clippy
            pkgs.rustfmt
            pkgs.rust-analyzer
            pkgs.cargo-nextest
            pkgs.pkg-config
            pkgs.git
          ];

          buildInputs = [ pkgs.zlib ];
        };

        benchmark = (pkgs.mkShell.override { stdenv = pkgs.gcc13Stdenv; }) {
          packages = benchmarkBins ++ [
            benchPython
            pkgs.pkg-config
            pkgs.cmake
            pkgs.gnumake
            pkgs.git
            pkgs.curl
            pkgs.fontconfig
            pkgs.poppler-utils
            pkgs.mupdf
          ];

          buildInputs = [
            pkgs.zlib
            pkgs.hdf5
          ];
        };

      };
    };
}
