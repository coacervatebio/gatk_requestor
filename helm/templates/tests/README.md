# Testing

Coacervate relies on a lot of moving parts, and will do so more in the future as the project matures. As such, testing of the default best-practices workflow focuses on an integration test for each step. This testing approach also leverages Helm's built-in chart testing command instead of re-inventing the wheel.

Each step-test follows how the coacervate-requestor is designed to function in the wild. A core set of dependencies and resources (e.g. the reference genome) are built into the requestor image. Customizing the workflow relies on mounting your own Snakefile, config, inputs and so on. 

Each test consists of a single Job object defining these different configurations, as it would when the requestor is run normally. The [test job](/helm/templates/tests/test_job.yaml), however, adds an init container and a test-checking container which setup inputs and compare outputs, respectively. The inputs and expected outputs for the default workflow are contained in a dedicated [test assets](/test_assets/) directory. These different assets are then "registered" in the [test values](/helm/test_values.yaml) yaml file along with any specific comparison command. This rather unorthodox approach means we can test the requestor the way it's meant to be used in production, with no mocking/stubbing, and a single template file. Not bad!

In order to run the test suite, you'll need to install the chart release with the additional test-values.yaml. Once that's all up and running, all that's needed is to run `helm test`.