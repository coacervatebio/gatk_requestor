from typing import List, Tuple
from flytekit.core.pod_template import PodTemplate
from kubernetes.client.models import V1PodSpec, V1Volume, V1Container, V1VolumeMount, V1EnvVar, V1PersistentVolumeClaimVolumeSource

yagna_requestor_ps = V1PodSpec(
    volumes=[
        V1Volume(
            name='yagna-storage',
            persistent_volume_claim=V1PersistentVolumeClaimVolumeSource(
                claim_name='yagna-pvc'
            )
        )
    ],
    containers=[
        V1Container(
            name='primary',
            image='docker.io/rwgrim/docker-noop',
            image_pull_policy='IfNotPresent',
            volume_mounts=[
                V1VolumeMount(
                    name='yagna-storage',
                    sub_path='yagna-socket',
                    mount_path='/yagna'
                ),
            ],
            env=[
                V1EnvVar(
                    name='YAGNA_API_URL',
                    value='http://yagna-service:7465'
                ),
                V1EnvVar(
                    name='GSB_URL',
                    value='unix:/yagna/yagna.sock'
                )
            ]
        )
    ]
)

yagna_requestor = PodTemplate(pod_spec=yagna_requestor_ps, labels={"pod_spec": "yagna_requestor"})
