# Using Singularity containers

Tested with `singularity-ce version 3.9.7-focal`

## Building containers

### KneadData


run

```sh
$ sudo singularity build $CONTAINER_PATH/kneaddata.sif $REPO_PATH/singularity_files/kneaddata.def
```

where `CONTAINER_PATH` is where you want the container to live,
and `REPO_PATH` is the path to this repository.

### MetaPhlAn

run

```sh
$ sudo singularity build $CONTAINER_PATH/metaphlan.sif $REPO_PATH/singularity_files/metaphlan.def
```

where `CONTAINER_PATH` is where you want the container to live,
and `REPO_PATH` is the path to this repository.

### HUMAnN

run

```sh
$ sudo singularity build $CONTAINER_PATH/humann.sif $REPO_PATH/singularity_files/humann.def
```

where `CONTAINER_PATH` is where you want the container to live,
and `REPO_PATH` is the path to this repository.


### PanPhlAn

run

```sh
$ sudo singularity build $CONTAINER_PATH/panphlan.sif $REPO_PATH/singularity_files/panphlan.def
```

where `CONTAINER_PATH` is where you want the container to live,
and `REPO_PATH` is the path to this repository.

## Running containers

To run a biobakery program, use `singularity exec`. For example, to run `humann`,

```sh
$ singularity exec $CONTAINER_PATH/humann.sif humann $ARGS...
```

I created a shell script and aliased it to `/usr/local/bin/humann`:

```sh
#!/bin/sh

singularity exec /murray/containers/humann.sif humann "$@"
```

 Mounting additional files systems

By default, singularity only mounts your working directory and its parents.
If you need to mount other file systems, use `--bind`. 

Eg. if you're running the container from your home folder,
but need to include `/some_filesystem`, run

```sh
$ singularity exec --bind /some_filesystem:/some_filesystem $CONTAINER_PATH/humann.sif humann $ARGS...
```
