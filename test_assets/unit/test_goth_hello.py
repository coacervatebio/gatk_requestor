import asyncio
from typing import AsyncIterable

from yapapi import Golem, Task, WorkContext
from yapapi.log import enable_default_logger
from yapapi.payload import vm

greeting = "\nHello from Goth!"


async def worker(context: WorkContext, tasks: AsyncIterable[Task]):
    async for task in tasks:
        script = context.new_script()
        future_result = script.run("/bin/echo", greeting)

        yield script

        task.accept_result(result=await future_result)


async def main():
    package = await vm.repo(
        image_hash="5e1013bb93f876d38b0fd6c24435bd118a5700c118a64c5d7f61630f",  # alpine vols only
    )

    tasks = [Task(data=None)]

    async with Golem(budget=1.0, subnet_tag="public") as golem:
        async for completed in golem.execute_tasks(worker, tasks, payload=package):
            print(completed.result.stdout)


def test_goth_hello():

    loop = asyncio.get_event_loop()
    task = loop.create_task(main())
    loop.run_until_complete(task)
