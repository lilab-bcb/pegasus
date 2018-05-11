FROM regevlab/scrtools:latest

RUN apt-get update && apt-get install --no-install-recommends -y \
    default-jdk

RUN mkdir /software/graph_layout/
ADD gephi-toolkit-0.9.2-all.jar /software/graph_layout/
ADD GraphLayout.class /software/graph_layout/

