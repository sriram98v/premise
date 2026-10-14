{
  #
  #     nix develop                    # dev shell: build and hack on PREMISE
  #     nix build                      # the PREMISE binary
  #     nix build .#premise-static     # Linux: a fully static (musl) binary that runs without Nix
  #     nix flake check                # build, tests, clippy and rustfmt, as CI runs them
  #
  # The benchmark against other classifiers, and its toolchain, live in the separate
  # premise-comparative-analysis repository, which takes this flake as an input.
  #
  description = "PREMISE — paired-end metagenomic read classifier";

  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});

      premiseFor = pkgs: pkgs.rustPlatform.buildRustPackage {
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
      packages = forAllSystems (pkgs:
        let premise = premiseFor pkgs; in {
          inherit premise;
          default = premise;
        } // pkgs.lib.optionalAttrs pkgs.stdenv.isLinux {
          premise-static = premiseFor pkgs.pkgsStatic;
        });

      checks = forAllSystems (pkgs:
        let premise = premiseFor pkgs; in {
          build = premise;

          tests = premise.overrideAttrs { doCheck = true; };

          # The benches are left out: they do not build against the current criterion.
          clippy = premise.overrideAttrs (old: {
            pname = "premise-clippy";
            nativeBuildInputs = old.nativeBuildInputs ++ [ pkgs.clippy ];
            buildPhase = ''
              runHook preBuild
              cargo clippy --release --lib --bins --tests -- -D warnings
              runHook postBuild
            '';
            installPhase = "touch $out";
          });

          fmt = pkgs.runCommand "premise-fmt" {
            nativeBuildInputs = [ pkgs.cargo pkgs.rustfmt ];
          } ''
            cd ${self}
            cargo fmt --check
            touch $out
          '';
        });

      devShells = forAllSystems (pkgs: {
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
      });
    };
}
