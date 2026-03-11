"""AWS Batch job submission utilities."""
import boto3
import json
import argparse
from pathlib import Path

BATCH_CLIENT = boto3.client('batch', region_name='us-west-2')
JOB_DEFS_DIR = Path(__file__).parent / 'job_definitions'


def submit_job(job_name, job_def_file, parameters=None, job_queue='perturb-seq-queue'):
    """Submit a job to AWS Batch.

    Parameters
    ----------
    job_name : str
    job_def_file : str  — filename in batch/job_definitions/
    parameters : dict, optional  — override job parameters
    job_queue : str

    Returns
    -------
    dict with jobId and jobArn
    """
    job_def_path = JOB_DEFS_DIR / job_def_file
    with open(job_def_path) as f:
        job_def = json.load(f)

    job_definition = job_def['jobDefinitionName']
    container_overrides = {}
    if parameters:
        container_overrides['environment'] = [
            {'name': k, 'value': str(v)} for k, v in parameters.items()
        ]

    response = BATCH_CLIENT.submit_job(
        jobName=job_name,
        jobQueue=job_queue,
        jobDefinition=job_definition,
        containerOverrides=container_overrides,
    )
    print(f"Submitted job: {response['jobId']}")
    return response


def main():
    parser = argparse.ArgumentParser(description='Submit AWS Batch jobs')
    parser.add_argument('job_type', choices=['scvi', 'de', 'cross_dataset'],
                        help='Job type to submit')
    parser.add_argument('--dataset', required=True, help='Dataset name or S3 path')
    parser.add_argument('--output-prefix', required=True, help='S3 output prefix')
    args = parser.parse_args()

    job_map = {
        'scvi': 'scvi_training.json',
        'de': 'genome_scale_de.json',
        'cross_dataset': 'cross_dataset_meta.json',
    }

    submit_job(
        job_name=f"{args.job_type}-{args.dataset.replace('/', '-')}",
        job_def_file=job_map[args.job_type],
        parameters={
            'DATASET': args.dataset,
            'OUTPUT_PREFIX': args.output_prefix,
        },
    )


if __name__ == '__main__':
    main()
