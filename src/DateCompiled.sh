#!/bin/sh
DATE=$(date)
DATE_SHORT=$(date +%Y-%m-%d)
VERSION=$(git log -1 --pretty="%h")
echo "const char *DateCompiled = \"$DATE (rev. $VERSION)\";" >DateCompiled.cpp
echo "const char *LogoDate = \"$DATE_SHORT (rev. $VERSION)\";" >>DateCompiled.cpp
