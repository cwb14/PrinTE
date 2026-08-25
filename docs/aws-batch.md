## Running on AWS Batch

The `awsbatch` profile runs the same pipeline on AWS with an S3 work directory. A sweep is
embarrassingly parallel, so Batch scales it about as well as anything will.

[Back to the README](../README.md)

---

> **I have not run this against a live AWS account.** The profile and the IAM policy below
> are what I would use, but nobody has paid a bill with them yet. Check your own costs
> before you launch a large sweep.

### What you need first

**An S3 bucket** for the work directory and the results. Nextflow stages inputs in and
outputs back out through it, so put it in the same region as the compute.

**A Batch compute environment.** Managed, `SPOT` if you can tolerate interruptions (the
pipeline retries), with an AMI that has the ECS agent and enough instance storage for a
genome plus its intermediates. `SIMULATE` is the process that matters: it wants the memory
described in [Troubleshooting](troubleshooting.md#how-much-memory-do-i-need), so choose an
instance family with at least 4 GB per vCPU. `r5` or `m5` are reasonable.

**A Batch job queue** attached to that compute environment. Its name is what you pass as
`--awsqueue`.

**An instance role** on the compute environment that can reach your bucket:

```json
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": ["s3:GetObject", "s3:PutObject", "s3:DeleteObject"],
      "Resource": "arn:aws:s3:::my-bucket/*"
    },
    {
      "Effect": "Allow",
      "Action": ["s3:ListBucket", "s3:GetBucketLocation"],
      "Resource": "arn:aws:s3:::my-bucket"
    }
  ]
}
```

The role also needs `AmazonECSTaskExecutionRolePolicy` so it can pull the container.

### Running it

```bash
nextflow run . -profile awsbatch,docker \
    -work-dir s3://my-bucket/printe-work \
    --awsqueue my-batch-queue \
    --awsregion us-east-1 \
    --mode sweep \
    --fasta s3://my-bucket/input/burnin.fasta \
    --bed   s3://my-bucket/input/burnin.bed \
    --clean_lib s3://my-bucket/input/lib_clean.fa \
    --ref_tsv   s3://my-bucket/input/real_LTR.tsv \
    --ref_fasta s3://my-bucket/input/real_genome.fa \
    --outdir s3://my-bucket/printe-results
```

`s3://` paths work anywhere a local path does, for inputs and for `--outdir`.

### Notes that will save you money or time

**Staging dominates, not compute.** A 3 Gb genome moving in and out of S3 for every task
costs more wall clock than the simulation. `nf/conf/awsbatch.config` sets
`maxParallelTransfers = 8`; raise it if your instances have the bandwidth, and prefer
fewer, longer generation steps (`--step`) over many short ones so each task does more work
per transfer.

**Use spot, and let retries do their job.** `nf/conf/base.config` retries on the exit codes
Batch uses for interruption and OOM, doubling memory each time.

**Clean the work directory.** `s3://my-bucket/printe-work` holds every intermediate genome
and will get large fast. `nextflow clean -f` after a successful run, or set an S3 lifecycle
rule to expire that prefix.

**Check the region.** Cross-region transfer charges are the usual surprise.

### Credentials

Nextflow picks up the standard chain: environment variables, `~/.aws/credentials`, or the
instance profile if you launch from EC2. Nothing PrinTE-specific. Do not put keys in
`nextflow.config`.
