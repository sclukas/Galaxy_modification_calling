# Galaxy Helm Modification Calling

Galaxy Docker image for RNA analysis and modification calling.

## Requirements
Make sure docker is installed on your computer.
```
sudo snap install docker
```

Access to the Galaxy Toolshed requires the addition of your nameserver to the /etc/resolv.conf file. This is necessary
for the first installation as some tools are downloaded from the Toolshed.
Get the nameserver with the following command:
```
nmcli dev show | grep DNS
```
Copy the nameserver and add it to the resolv.conf file. You can find the file in the /etc directory.
Add "nameserver <xx.xx.xxx.x>" and save the file (this might require root rights).
**The resolv.conf file resets after every restart. If you have to install programs from the Toolshed you have to add the nameserver again**

## Setup
Enter your Galaxy-folder and open a terminal within the directory.
Use the following command to create the Docker image:
```
sudo docker build -t galaxy-helm setup
```

## Run Galaxy
A fresh instance of the created image can now be started using the startup script:
```
sudo ./run.sh
```
You may have to make the file (run.sh) executable.

The docker image and the Galaxy instance require about 11 gigabytes of memory.
After a short startup period, the Galaxy instance should be accessible via <http://localhost:8080>.

## First run
The Galaxy instance provides a single admin user with user name `admin@galaxy.org` and password `admin` that can be used to create additional user accounts.

**You might want to consider changing the default password as soon as possible!**

That's it! You can now upload your data and create workflows (or import existing ones from the `workflow/` directory).


## Additional useful commands


# List all containers
sudo docker ps -a
# Delete container
sudo docker rm CONTAINER_ID
# Stop container
docker stop CONTAINER_NAME

# List images
sudo docker image ls
# Remove Docker image
sudo docker rmi IMAGE_ID
sudo docker info

# Restart docker (free port)
sudo service docker restart

# restart Galaxy without restarting the container
sudo docker exec <container name> supervisorctl restart galaxy:

# Start shell in docker
sudo docker exec -it "id of running container" bash
