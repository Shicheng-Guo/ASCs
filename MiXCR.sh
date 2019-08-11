# install ruby

sudo apt update
sudo apt install ruby-full
ruby --version

# install brew
sh -c "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install.sh)"

# install JDK
# GO TO: https://www.oracle.com/technetwork/java/javase/downloads/jdk12-downloads-5295953.html
which java
whereis java
vi ~/.bash_profile
export JAVA_HOME=/home/guosa/hpc/tools/jdk-12.0.2/
source ~/.bash_profile
echo $JAVA_HOME
which java

# install maven
wget http://apache.mirrors.tds.net/maven/maven-3/3.6.1/binaries/apache-maven-3.6.1-bin.tar.gz
export PATH=/home/guosa/hpc/tools/apache-maven-3.6.1/bin:$PATH
# Error: Could not find or load main class org.codehaus.plexus.classworlds.launcher.Launcher
# If you meet above problem which indicates you downloaded wrong file apache-maven-3.6.1-src.tar.gz. 
# What you need to do is just re-download  apache-maven-3.6.1-bin.tar.gz 

# install mixcr
wget https://github.com/milaboratory/mixcr/releases/download/v3.0.9/mixcr-3.0.9.zip
unzip mixcr-3.0.9.zip
