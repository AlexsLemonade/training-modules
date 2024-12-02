The shell script in this directory downloads processed `SingleCellExperiment` object for the `SCPCS000203` sample from `SCPCP000006` (Wilms tumor) using the data download mechanism from OpenScPCA.

You may want to run it with your OpenScPCA conda environment activated.

By default, it will use your currently active AWS profile (falling back to one called `openscpca` if not set) and download data from the `2024-11-25` OpenScPCA release.

You can alter the AWS profile or release with the following:

```sh
AWS_PROFILE={profile} RELEASE={release} ./download-openscpca-data.sh
```

Replacing `{profile}` and `{release}` with a profile with OpenScPCA access and valid release, respectively.
