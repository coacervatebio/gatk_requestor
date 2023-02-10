import pytest
import kubernetes
from kubernetes import client

@pytest.fixture
def k8s_api_client():
    # Load kubernetes configuration
    kubernetes.config.load_kube_config()
    # Create a Kubernetes API client
    api_client = client.BatchV1Api()
    return api_client

def test_run_kubernetes_jobs(k8s_api_client):
    # Define a list of commands to run as part of the Kubernetes jobs
    commands = [
        "echo 'Job 1'",
        "echo 'Job 2'",
        "echo 'Job 3'"
    ]

    # Loop through the list of commands and trigger a Kubernetes job for each command
    for command in commands:
        job = client.V1Job(
            api_version="batch/v1",
            kind="Job",
            metadata=client.V1ObjectMeta(name="test-job"),
            spec=client.V1JobSpec(
                template=client.V1PodTemplateSpec(
                    metadata=client.V1ObjectMeta(labels={"app": "test-job"}),
                    spec=client.V1PodSpec(
                        containers=[
                            client.V1Container(
                                name="test-container",
                                image="alpine",
                                command=["/bin/sh"],
                                args=["-c", command],
                            )
                        ],
                        restart_policy="Never"
                    )
                )
            )
        )
        # Create the Kubernetes job
        k8s_api_client.create_namespaced_job(namespace="default", body=job)
        # Wait for the job to complete
        k8s_api_client.read_namespaced_job_status(name="test-job", namespace="default")
        # Delete the job once it has completed
        k8s_api_client.delete_namespaced_job(name="test-job", namespace="default", body=client.V1DeleteOptions(propagation_policy='Foreground', grace_period_seconds=5))
