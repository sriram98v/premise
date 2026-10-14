{
  #
  #     nix develop                    # dev shell: build and hack on PREMISE
  #     nix build                      # the PREMISE binary
  #
  # The benchmark against other classifiers, and its toolchain, live in the separate
  # premise-comparative-analysis repository, which takes this flake as an input.
  #
  description = "PREMISE — paired-end metagenomic read classifier";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      system = "x86_64-linux";
      pkgs = import nixpkgs { inherit system; };

      premise-pkg = pkgs.rustPlatform.buildRustPackage {
        pname = "premise";
        version = "0.5.0";
        src = self;

        cargoLock = {
          lockFile = ./Cargo.lock;
          outputHashes = {
            "haystackfm-0.6.8" = "sha256-1knSJSCm1SK0LTkH3JiCMyZ5Jo25rjBGxl7YR26ryYQ=";
          };
        };

        buildType = "release";
        doCheck = false;

        meta = with pkgs.lib; {
          description = "PREMISE: A probabilistic framework for source assignment of Illumina reads";
          homepage = "https://github.com/sriram98v/premise";
          license = licenses.mit;
          mainProgram = "premise";
          platforms = platforms.unix;
        };
      };
    in {
      packages.${system} = {
        premise = premise-pkg;
        default = premise-pkg;
      };

      devShells.${system}.default = pkgs.mkShell {
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
    };
}
