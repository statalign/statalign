#!/bin/bash -eu

# Don't forget to source your ~/.bashrc file afterwards!

set -o pipefail
this_dir="$(dirname "$0")"

gradle_version="2.0"
gradle_download_url="https://services.gradle.org/distributions/gradle-${gradle_version}-all.zip"
gradle_install_dir="$HOME/gradle-${gradle_version}-install"
gradle_home="${gradle_install_dir}/gradle-${gradle_version}/"


escape_string () 
{
    string="$1";
    echo -n "'$(echo -n "$string" | sed "s|'|'\\\''|g")'"
}

function check_command {
  declare command_name
  for command_name in "$@"; do
    compgen -c "$command_name" | grep -Fxqe "$command_name" || (echo "$command_name: command not found, but needed"; false )
  done
}

check_command "unzip" "java"

mkdir -p "${gradle_install_dir}"
pushd "${gradle_install_dir}" > "/dev/null"
  if [ ! -f ${gradle_install_dir}/gradle-${gradle_version}-all.zip ]; then
    echo "Downloading gradle ${gradle_version} ..."
    wget "${gradle_download_url}"
  fi
  echo "Extracting package..."
  unzip -o "$(basename "${gradle_download_url}")" > "/dev/null"
  echo ...done
popd > "/dev/null"

mkdir -p ${gradle_home}

if ! [ -d "$gradle_home" ]; then (
  echo "Installation not complete: $gradle_home is not a valid directory. PATH variable won't be set"
  false
) fi

# maybe it does not exist yet
touch ~/.bashrc

d_bashrc_lines="$(cat << BASHRC_LINES
export GRADLE_HOME=$(escape_string "$gradle_home")
export PATH=$(escape_string "$gradle_home/bin"):"\$PATH"
BASHRC_LINES
)"

eval "$d_bashrc_lines"
gradle -version
echo "$d_bashrc_lines" >> "$HOME/.bashrc_gradle"

bashrc_line='. "$HOME/.bashrc_gradle" # automatic generated line w3'
if ! grep -Fxqe "$bashrc_line" "$HOME/.bashrc"; then
  cat "$HOME/.bashrc" > "$HOME/.bashrc_tmp"
  echo "$bashrc_line" >> "$HOME/.bashrc_tmp"
  mv "$HOME/.bashrc_tmp" "$HOME/.bashrc" 
fi
