FROM python:3.12-slim
WORKDIR /app
COPY server.py ./
COPY docs ./docs
ENV HOST=0.0.0.0 PORT=8000
EXPOSE 8000
USER nobody
CMD ["python3", "server.py"]
