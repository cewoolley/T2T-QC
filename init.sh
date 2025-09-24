# init docker etc

sudo yum install -y yum-utils
sudo yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
sudo yum install -y docker-ce docker-ce-cli containerd.io

sudo systemctl enable docker
sudo systemctl start docker

docer --version
sudo systemctl status docker

docker run hello-world

sudo groupadd docker
sudo usermod -aG docker datascience
newgrp docker

docker run hello-world

