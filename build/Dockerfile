FROM python:3.9-slim

# Set the working directory in the container
WORKDIR /proj/gibbons/nbohmann/personalized_predictions/build

# Copy requirements.txt and install Python dependencies
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

# Copy the CPLEX distribution into the Docker image
COPY /proj/gibbons/software/ibm/cplex/python/3.8/x86-64_linux /cplex

# Install the CPLEX distribution
RUN pip install /cplex

# Copy the current directory contents into the container at /usr/src/app
COPY . .

# Make port 80 available to the world outside this container
EXPOSE 80

# Run app.py when the container launches
CMD ["python", "build.py"]