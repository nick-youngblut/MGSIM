# Use an official Python runtime as a base image
FROM ubuntu:20.04

# Set the environment variable to non-interactive
ENV DEBIAN_FRONTEND noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    build-essential \
    curl \
    vim \
    git \
    less \
    tree \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Reset the environment variable
ENV DEBIAN_FRONTEND=

