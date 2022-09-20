FROM python:3.10.7-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest

WORKDIR /app

COPY requirements.txt /app/requirements.txt

RUN pip install -r requirements.txt

COPY ./R/requirements.r /app/requirements.r

RUN Rscript requirements.r

# COPY . /app