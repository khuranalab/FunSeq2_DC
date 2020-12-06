### Load and Run Funseq2 Docker Container Instructions

#### 1. Download the docker desktop:
- Visit the docker official site at https://www.docker.com/get-started
- Download the docker desktop for your operating system

#### 2. Download the data context files:
- Funseq2 requires data context files which can be obtained from https://khuranalab.med.cornell.edu/data_DC3.html


#### 3. Load the docker image:

- Load the docker image from the tar archive using the following command:

    - For macOS and Linux:
    
        ```bash
        docker load < funseq2-docker-latest.tar
        ```   
    - For Windows:
        ```bash
        docker load -i funseq2-docker-latest.tar
        ``` 
- NOTE: Depending on your system privilege level, you may need to prefix this command with ```sudo```


#### 4. Modify the paths:
- Open ```funseq2-docker-latest.sh``` script and modify the paths for the ```data_context``` and ```user_input``` to your local paths before running the script
- Run the following command:

    ```bash
    ./funseq2-docker-latest.sh
    ```
