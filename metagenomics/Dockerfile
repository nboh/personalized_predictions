FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /proj/gibbons/nbohmann/personalized_predictions/metagenomics

# Copy the conda environment file into the container
COPY conda.yml ./

# Create the conda environment
RUN conda env create -f conda.yml

# Activate the environment
SHELL ["conda", "run", "-n", "metagenomics", "/bin/bash", "-c"]

# Copy the current directory contents into the container
COPY . .

# Make port 80 available to the world outside this container
EXPOSE 80

# Install Nextflow
RUN apt-get update && \
    apt-get install -y openjdk-11-jdk wget && \
    wget -qO- https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/

# Run the Nextflow pipeline when the container launches
CMD ["nextflow", "run", "main.nf", "--data_dir=./data", "--single_end=false", "--refs=/proj/gibbons/refs", "--trim_front=5", "--min_length=50", "--quality_threshold=20", "--read_length=150", "--threshold=10"]
