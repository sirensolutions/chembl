FROM ubuntu:latest
MAINTAINER Andrea Pierleoni "andrea.pierleoni@siren.io"
RUN apt-get update -y
RUN apt-get install -y python-pip python-dev build-essential python-rdkit librdkit1 rdkit-data
COPY . /app
WORKDIR /app
RUN pip install -r requirements.txt
#ENTRYPOINT ["python"]
EXPOSE 8080
CMD ["gunicorn", "-b", ":8080", "fingerprint_api:app"]
