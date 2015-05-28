#include "fs_tcp_client.h"

#include <iostream>
#include <QSharedPointer>

namespace fs {

#define TCPIPMAGICNUMBER 33433

/////////////////////////////////////////////////////////////////////////////////////
FSTCPClient::FSTCPClient (QObject *parent)
  : QObject (parent)
  , tcp_socket_ (0)
{

}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::init ()
{
  std::cerr << "FSTCPClient init...\n";
  tcp_socket_ = new QTcpSocket();
  connect (tcp_socket_, SIGNAL(error(QAbstractSocket::SocketError)), SLOT(error(QAbstractSocket::SocketError)));
  connect (tcp_socket_, SIGNAL(disconnected()), SLOT(disconnected()) );
  connect (tcp_socket_, SIGNAL(connected()), SLOT(connected()) );
  connect (tcp_socket_, SIGNAL(readyRead()), SLOT(readyRead()));
  connect (tcp_socket_, SIGNAL(stateChanged(QAbstractSocket::SocketState)), SLOT(stateChanged(QAbstractSocket::SocketState)) );

  message_block_size_ = 0;
  callback_object_ = 0;
  new_message_ = true;

  decoder_ = new fsBinaryStream ();
}


/////////////////////////////////////////////////////////////////////////////////////
FSTCPClient::~FSTCPClient ()
{
}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::setCallbackObject (FaceShiftClient *obj)
{
  callback_object_ = obj;
}


/////////////////////////////////////////////////////////////////////////////////////
bool
FSTCPClient::isConnected ()
{
  if (tcp_socket_->state () == QAbstractSocket::ConnectedState)
    return true;
  return false;
}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::readyRead ()
{
  if (tcp_socket_->state () != QAbstractSocket::ConnectedState)
    return;

  // data to read on listen socket
  if (tcp_socket_->isValid ())
  {
    char data[2048];
    qint64 length = tcp_socket_->read (data, 2048);
    decoder_->received (length, data);

    fsMsgPtr message;
    while ((message = decoder_->get_message ()))
    {
      /// callback call should go here.
      if (message->id () == fsMsg::MSG_OUT_TRACKING_STATE)
      {
        QSharedPointer<fsMsgTrackingState> new_state (new fsMsgTrackingState (*dynamic_cast<fsMsgTrackingState*> (message.data ())));
        emit receivedNetworkUpdates (new_state);
      }
    }
  }
  else
  {
    if (tcp_socket_->state () != QAbstractSocket::ConnectedState)
      disconnect ();
  }
}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::error (QAbstractSocket::SocketError error)
{
  disconnect();
  if (callback_object_)
    callback_object_->setNetworkConnectionState (false);
}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::stateChanged (QAbstractSocket::SocketState state)
{

}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::connected ()
{
  if (callback_object_)
    callback_object_->setNetworkConnectionState (true);
}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::disconnected()
{
  if (tcp_socket_->state () == QAbstractSocket::ConnectedState)
  {
    message_block_size_ = 0;
    new_message_ = true;
  }

  if (callback_object_)
    callback_object_->setNetworkConnectionState (false);
}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::connectToServer (QString hostname, int port)
{
  while (!tcp_socket_);
  printf ("connecting to server %s on port %d...\n", hostname.toStdString ().c_str (), port);
  disconnect();

  message_block_size_ = 0;
  new_message_ = true;
  tcp_socket_->connectToHost (hostname, port);
  if (!tcp_socket_->waitForConnected (1000))
  {
    disconnect ();
    callback_object_->setNetworkConnectionState (false);
  }
}


/////////////////////////////////////////////////////////////////////////////////////
void
FSTCPClient::disconnect ()
{
  if (tcp_socket_->state () == QAbstractSocket::ConnectedState)
  {
    tcp_socket_->close ();
    message_block_size_ = 0;
    new_message_ = true;
  }
}


}
