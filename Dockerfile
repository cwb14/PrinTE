# Build stage: solve the environment and compile the point-mutation binary.
FROM mambaorg/micromamba:1.5-jammy AS build

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ make \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -afy

COPY --chown=$MAMBA_USER:$MAMBA_USER . /src
WORKDIR /src
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN make ltr-mutator PREFIX=/tmp/build \
    && pip install --no-deps --no-cache-dir .

# Runtime stage.
FROM mambaorg/micromamba:1.5-jammy

LABEL org.opencontainers.image.title="PrinTE" \
      org.opencontainers.image.description="Forward simulator of transposable-element genome evolution" \
      org.opencontainers.image.source="https://github.com/cwb14/PrinTE" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

COPY --from=build /opt/conda /opt/conda
COPY --from=build /tmp/build/bin/ltr_mutator /opt/conda/bin/ltr_mutator
COPY --from=build /src/PrinTE.sh /opt/printe/PrinTE.sh
COPY --from=build /src/data /opt/printe/data

# The image layers are read-only at runtime, so the mutator is already on PATH and
# the Kmer2LTR clone that post-processing performs goes to a writable temp dir.
ENV PRINTE_MUTATOR=/opt/conda/bin/ltr_mutator \
    PRINTE_SCRIPT=/opt/printe/PrinTE.sh \
    PRINTE_DATA=/opt/printe/data \
    PRINTE_CACHE=/tmp/printe \
    PATH=/opt/conda/bin:$PATH

USER $MAMBA_USER
WORKDIR /work
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["printe"]
