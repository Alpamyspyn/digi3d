#include "fsbinarystream.h"
#include <fstream>

bool
loadFaceshiftAnimation (std::string const &path,
                        std::vector<fs::fsMsgTrackingState> &tracking_data)
{
  std::ifstream stream;
  stream.open (path.c_str (), std::ios_base::in | std::ios_base::binary);

  if (!stream.is_open ())
    return (false);

  long int const bufferSize = 512;
  long int read;
  char buffer[bufferSize];
  fs::fsBinaryStream fsstream;
  do {
    stream.read (buffer, bufferSize);
    read = stream.gcount ();
    // Send this segment over to Faceshift for decoding
    fsstream.received (read, buffer);

    fs::fsMsgPtr message;
    do {
      message = fsstream.get_message();
      if(message != NULL && message->id() == fs::fsMsg::MSG_OUT_TRACKING_STATE)
        tracking_data.push_back (* dynamic_cast<fs::fsMsgTrackingState *>(message.data ()) );
    } while(fsstream.valid() && message != NULL);

  } while (stream.good ());


  stream.close ();
  return (true);
}

