{
  # premise, plus the pinned toolchain for its comparative benchmark.
  #
  #     nix develop                    # dev shell: build and hack on premise
  #     nix develop .#benchmark        # study shell: every tool comparative-analysis/run-analysis.py needs
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
        version = "0.4.0-${self.shortRev or self.dirtyShortRev or "dirty"}";
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
        src = pkgs.fetchurl {
          url = "https://github.com/bwa-mem2/bwa-mem2/releases/download/v${version}/Source_code_including_submodules.tar.gz";
          hash = "sha256-DEih6oAK9JmucmS0yJCMTKNfvlp98q1hBseaqUu0nLs=";
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

      cnsplotsPkgs = pkgs.callPackage ./nix/cnsplots.nix { };

      benchPython = pkgs.python3.withPackages (ps: [
        ps.matplotlib
        ps.pandas
        ps.pysam
        fastapy
        (py.toPythonModule insilicoseq)
        cnsplotsPkgs.cnsplots
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

        bwa-mem2
        cutadapt
        insilicoseq
        pkgs.sratoolkit
        pkgs.samtools
        pkgs.seqkit

        pkgs.corefonts
      ];
    in {
      packages.${system} = {
        inherit cutadapt dnaio xopen fastapy bwa-mem2 insilicoseq;
        inherit (methods) raptor ganon ganon-cpp multitax pylowestcommonancestor karp centrifuger mora sylph kmcp;
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

          shellHook = ''
            echo "premise dev shell"
            echo "  rustc     $(rustc --version 2>/dev/null | awk '{print $2}')"
            echo "  cargo     $(cargo --version 2>/dev/null | awk '{print $2}')"
            echo ""
            echo "  cargo build --release   |   cargo test   |   cargo bench"
            echo "  benchmark toolchain:  nix develop .#benchmark"
          '';
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

          shellHook = ''
            export CMAKE_PREFIX_PATH="${pkgs.zlib.dev}:${pkgs.zlib}''${CMAKE_PREFIX_PATH:+:$CMAKE_PREFIX_PATH}"
            export LC_ALL=C

            export FONTCONFIG_FILE=${benchFontsConf}
            # matplotlib invalidates its font cache on a matplotlib version bump and on nothing
            # else -- never on newly installed fonts -- so a stale cache would hide Arial and trip
            # ablation_figures.py's guard. Keying the cache dir on the fonts.conf store path makes
            # it self-bust whenever the font set changes.
            export MPLCONFIGDIR="''${XDG_CACHE_HOME:-$HOME/.cache}/premise-mpl/${builtins.baseNameOf benchFontsConf}"
            mkdir -p "$MPLCONFIGDIR"

            echo "premise benchmark shell   (every binary below comes from this flake)"
            echo "  cc        $(cc --version 2>/dev/null | head -1)"
            echo "  premise   ${premise-pkg.version}"
            echo "  methods   karp centrifuger mora sylph kmcp ganon-${methods.ganon.version} raptor-${methods.raptor.version}"
            echo "  aux       bwa-mem2-${bwa-mem2.version} samtools-${pkgs.samtools.version} cutadapt-${cutadapt.version} seqkit-${pkgs.seqkit.version}"
            echo "  data      iss-${insilicoseq.version} sratoolkit-${pkgs.sratoolkit.version}"
            echo "  python3   $(python3 --version 2>&1 | awk '{print $2}')  (stdlib + matplotlib, pandas, pysam, fastapy, iss)"
            echo "  fonts     Arial (corefonts) via FONTCONFIG_FILE; figure audit: pdffonts, mutool"
          '';
        };

      };
    };
}
