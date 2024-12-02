The shell script in this directory downloads processed `SingleCellExperiment` objects and metadata from tumor samples in `SCPCP000015` using the data download mechanism from OpenScPCA.

You may want to run it with your OpenScPCA conda environment activated.

By default, it will use an AWS profile called `openscpca` and download data from the `2024-11-25` OpenScPCA release.

You can alter the AWS profile or release with the following:

```sh
PROFILE={profile} RELEASE={release} ./download-openscpca-data.sh
```

Replacing `{profile}` and `{release}` with a profile with OpenScPCA access and valid release, respectively.
