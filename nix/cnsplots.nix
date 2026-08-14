# cnsplots — Cell/Nature/Science figure styling, used by comparative-analysis/scripts/ablation_figures.py.
#
{ lib, python3Packages, fetchPypi }:

let
  ps = python3Packages;

  simple = { pname, version, hash, deps ? [ ] }:
    ps.buildPythonPackage {
      inherit pname version;
      pyproject = true;
      src = fetchPypi { inherit pname version hash; };
      build-system = [ ps.setuptools ];
      dependencies = deps;
      doCheck = false;
      pythonImportsCheck = [ pname ];
    };

in rec {
  num2tex = simple {
    pname = "num2tex";
    version = "0.8";
    hash = "sha256-1bXirveaFBk0yCTTiX+d2f250l58BCu8AWHnZHzWb5k=";
  };

  statannotations = simple {
    pname = "statannotations";
    version = "0.7.2";
    hash = "sha256-qT5ChoLWbJ/KWjmdW9P+fJzvBQD1OEDG9TYlJrRQHpM=";
    deps = [ ps.matplotlib ps.numpy ps.pandas ps.scipy ps.seaborn ];
  };

  cnsplots = ps.buildPythonPackage rec {
    pname = "cnsplots";
    version = "0.5.0";
    pyproject = true;
    src = fetchPypi {
      inherit pname version;
      hash = "sha256-oCuXSsDPD/xGykf1Z9qAeSBznhsN4QoItnzcU82ms2o=";
    };
    build-system = [ ps.setuptools ];

    dependencies = [
      ps.lxml
      ps.matplotlib
      ps.numpy
      ps.palettable
      ps.pandas
      ps.scipy
      ps.seaborn
      num2tex
      statannotations
    ];

    pythonRemoveDeps = [
      "adjustText"
      "anndata"
      "biopython"
      "comprisk"
      "gseapy"
      "lifelines"
      "matplotlib-venn"
      "natsort"
      "patsy"
      "pycomplexheatmap"
      "scanpy"
      "scikit-learn"
      "statsmodels"
      "upsetplot"
    ];

    pythonRelaxDeps = [ "matplotlib" ];

    doCheck = false;
    pythonImportsCheck = [ "cnsplots" ];

    meta = with lib; {
      description = "Publication-style plotting helpers (Cell/Nature/Science figure conventions)";
      homepage = "https://pypi.org/project/cnsplots/";
      license = licenses.mit;
    };
  };
}
