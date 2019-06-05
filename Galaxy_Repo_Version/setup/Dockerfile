# Galaxy - Helm pipeline
#
# VERSION 1.0

FROM bgruening/galaxy-stable:19.01

WORKDIR /galaxy-central


# Additional shed tools
RUN pip install --upgrade pip

ADD shed-tools-extra.yml $GALAXY_ROOT/shed-tools-extra.yml
RUN install-tools $GALAXY_ROOT/shed-tools-extra.yml

# Local tools
ADD local-tools.xml /local-tools/tools.xml
ADD tools/* /local-tools/tools/
ENV GALAXY_CONFIG_TOOL_CONFIG_FILE config/tool_conf.xml.main,config/shed_tool_conf.xml,/local-tools/tools.xml

RUN pip install numpy scipy pandas matplotlib scikit-learn
#RUN rm -f /usr/bin/python && ln -s /usr/bin/python /usr/bin/python3
ADD extra_files/__init__.py $GALAXY_ROOT/lib/galaxy/util/__init__.py
ADD extra_files/binary.py $GALAXY_ROOT/lib/galaxy/datatypes/binary.py
ADD extra_files/registry.py $GALAXY_ROOT/lib/galaxy/datatypes/registry.py
ADD extra_files/datatypes_conf.xml $GALAXY_ROOT/config/datatypes_conf.xml

# Customize some Galaxy settings
ENV GALAXY_CONFIG_ALLOW_USER_CREATION True
ENV GALAXY_CONFIG_ALLOW_USER_DELETION True
ENV GALAXY_CONFIG_BRAND "Helm Pipeline"
ENV GALAXY_CONFIG_REQUIRE_LOGIN True

# Add workflows to the Docker image
#ADD ./workflows/*.ga $GALAXY_ROOT/workflows/


# Mark folders as imported from the host.
VOLUME ["/export/", "/data/", "/var/lib/docker"]


# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :443
EXPOSE :8800

# Autostart script that is invoked during container start
CMD ["/usr/bin/startup"]
