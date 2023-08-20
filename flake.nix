{
  description = "Description for the project";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    devenv.url = "github:cachix/devenv";
    # nix2container.url = "github:nlewo/nix2container";
    # nix2container.inputs.nixpkgs.follows = "nixpkgs";
    mk-shell-bin.url = "github:rrbutani/nix-mk-shell-bin";
    dream2nix.url = "github:nix-community/dream2nix";
  };

  outputs = inputs @ {
    flake-parts,
    dream2nix,
    devenv,
    ...
  }:
    flake-parts.lib.mkFlake {inherit inputs;} {
      imports = [
        inputs.devenv.flakeModule
        inputs.dream2nix.flakeModuleBeta
      ];
      systems = ["x86_64-linux"];

      perSystem = {
        config,
        self',
        inputs',
        pkgs,
        system,
        ...
      }: {
        formatter = pkgs.alejandra;

        devenv.shells.default = {
          name = "selfingsim";

          # languages.python.enable = true;

          # https://devenv.sh/reference/options/
          packages = with pkgs; [
            hatch 
            (python310.withPackages (ps: with ps; [hatchling]))
            zlib
            gcc
            nodejs-slim
            ripgrep
          ];

          pre-commit.hooks = {
            black.enable = true;
            isort.enable = true;
          };

          enterShell = ''
            echo "Hello"
          '';
        };

        dream2nix.inputs."selfingsim" = {
          source = ./.;
          projects.selfingsim = {
            subsystem = "python";
            translator = "pip";
            subsystemInfo.system = system;
            subsystemInfo.pythonVersion = "3.10";
          };
          packageOverrides = {
            selfingsim = {
              "simupop-deps" = {
                overrideAttrs = oldAttrs: {
                  buildInputs = oldAttrs.buildInputs ++ [pkgs.zlib.dev];
                };
              };
            };
          };
        };

        packages = {
          inherit (config.dream2nix.outputs."selfingsim".packages) selfingsim resolveImpure;
          inherit (inputs'.devenv.packages) devenv;
        };
      };
      flake = {
        # The usual flake attributes can be defined here, including system-
        # agnostic ones like nixosModule and system-enumerating ones, although
        # those are more easily expressed in perSystem.
      };
    };
}
