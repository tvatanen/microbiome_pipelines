# The build-stage image:
FROM mambaorg/micromamba:1.4.3

# cd condaenvs/humann
# sudo docker build -t humann-local .
# docker run -it --rm humann-local humann --h

# sudo docker build -t humann-nodb:mamba-v3.7 .
# docker tag humann-nodb:mamba-v3.7 public.ecr.aws/j5i5h1i5/humann-nodb:mamba-v3.7
# docker push public.ecr.aws/j5i5h1i5/humann-nodb:mamba-v3.7

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
