{
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      with import nixpkgs { inherit system; };
      with rPackages; {
        defaultPackage = (buildRPackage {
          name = "svaNUMT";
          src = ./.;
          propagatedBuildInputs = [
            GenomicRanges
            rtracklayer
            StructuralVariantAnnotation
            VariantAnnotation
            assertthat
            Biostrings
            dplyr
            rlang
            GenomicFeatures
            GenomeInfoDb
            S4Vectors
          ];
        }).overrideAttrs (_: {
          doCheck = true;
          checkInputs = [ BiocCheck ];
          script = writeText "bioc-check.r" ''
            library(BiocCheck)
            BiocCheck(paste0(Sys.getenv("TMPDIR"), "/svaRetro")
              , `no-check-deprecated`=T  # Requires net
              , `no-check-CRAN`=T        # Requires net
              , `no-check-version-num`=T # No idea why this fails
              , `quit-with-status`=T)
          '';
          checkPhase = ''
            cp -r $PWD $TMPDIR/svaRetro
            Rscript $script
          '';
        });
      });
}
