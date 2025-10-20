# SLURM Job Submission Guide

## Basic SLURM Commands

### Submit a job
```bash
sbatch script.sh
```

### Check job status
```bash
squeue -u $USER
```

### Cancel a job
```bash
scancel <job_id>
```

### View job details
```bash
scontrol show job <job_id>
```

### Check past jobs
```bash
sacct -u $USER
```

## Resource Allocation Guidelines

### For Pore-C Analysis

**Small datasets (<5GB fastq):**
```bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=6:00:00
```

**Medium datasets (5-20GB fastq):**
```bash
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=12:00:00
```

**Large datasets (>20GB fastq):**
```bash
#SBATCH --cpus-per-task=76
#SBATCH --mem=256G
#SBATCH --time=24:00:00
```

### For Hi-C Analysis

**Standard paired-end Hi-C:**
```bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=8:00:00
```

## Example SLURM Script

See `examples/slurm_submission_example.slurm` for a complete working example.

## Monitoring Jobs

### Check memory usage
```bash
sstat -j <job_id> --format=MaxRSS,AvgRSS
```

### Check output in real-time
```bash
tail -f slurm-<job_id>.out
```

## Common Issues

**Out of Memory:**
- Increase `--mem` allocation
- For Pore-C: reduce `--chunk_size` in nextflow command

**Timeout:**
- Increase `--time` allocation
- Consider splitting large jobs

**Job Pending:**
- Check queue: `squeue`
- Check account allocation: `mybalance` or equivalent
