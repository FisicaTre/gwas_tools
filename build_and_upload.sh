upload_type=$1

rm -rf build dist gwas_tools.egg-info
python3 -m build --wheel

if [[ "$upload_type" = "test" ]]; then
  python3 -m twine upload --verbose --skip-existing dist/* -r testpypi
elif [[ "$upload_type" = "prod" ]]; then
  python3 -m twine upload --verbose --skip-existing dist/*
fi

exit 0
