# Build stage: solve the environment and compile the point-mutation binary.
FROM mambaorg/micromamba:1.5-jammy AS build

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ make git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean -afy

COPY --chown=$MAMBA_USER:$MAMBA_USER . /src
WORKDIR /src
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN make ltr-mutator PREFIX=/tmp/build \
    && pip install --no-deps --no-cache-dir .

# Build in the Kmer2LTR commit that PrinTE.sh pins, so post-processing never clones it
# at runtime and works without network access.
RUN ref="$(sed -n 's/^KMER2LTR_REF=//p' PrinTE.sh)" \
    && test -n "$ref" \
    && git clone --quiet https://github.com/cwb14/Kmer2LTR.git /tmp/build/Kmer2LTR \
    && git -C /tmp/build/Kmer2LTR checkout --quiet "$ref" \
    && rm -rf /tmp/build/Kmer2LTR/.git

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
COPY --from=build /src/R /opt/printe/R
COPY --from=build /tmp/build/Kmer2LTR /opt/printe/Kmer2LTR

# The image layers are read-only at runtime, so everything PrinTE would otherwise build
# or fetch is already here: the mutator on PATH, Kmer2LTR next to PrinTE.sh. Setting
# PRINTE_CACHE to that same directory keeps a value from the host (Apptainer passes the
# host environment in) from sending PrinTE off to clone Kmer2LTR again.
ENV PRINTE_MUTATOR=/opt/conda/bin/ltr_mutator \
    PRINTE_SCRIPT=/opt/printe/PrinTE.sh \
    PRINTE_DATA=/opt/printe/data \
    PRINTE_CACHE=/opt/printe \
    PATH=/opt/conda/bin:$PATH

USER $MAMBA_USER
WORKDIR /work
# printe belongs in the ENTRYPOINT, not the CMD: arguments replace CMD, so with
# CMD ["printe"] a `docker run printe --version` would drop printe and try to exec
# --version. Putting it here makes the image runnable as a command, which is also
# what lets `./printe.sif --burnin_only ...` work under Apptainer.
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "printe"]
