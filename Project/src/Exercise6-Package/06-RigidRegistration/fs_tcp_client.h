#pragma once

#include <QObject>
#if (QT_MAJOR_VERSION == 4)
  #include <QTCPSocket>
#else
  #include <QtNetwork/QTcpSocket>
#endif

#include "fsbinarystream.h"


class FaceShiftClient : public QObject
{
Q_OBJECT
public:
  FaceShiftClient () {}

  void
  setNetworkConnectionState (bool) {}

public slots:
    void
    initializeConnection () {}

    void
    receivedNetworkUpdates (QSharedPointer<fs::fsMsgTrackingState>) {};

signals:
    void
    askForConnection (QString, int);
};

namespace fs {

/** \brief TCP client that communicates with the FaceShift Studio application. */
class FSTCPClient : public QObject
{
  Q_OBJECT
public:
  /** \brief Constructor that initializes the socket and the signals received from it. */
  explicit FSTCPClient (QObject *parent = 0);

  /** \brief Empty destructor. */
  virtual ~FSTCPClient ();

  /** \brief Disconnect current connection. */
  void
  disconnect ();

  /** \brief Method to set the object to which the callbacks will be made.
      * \param[in] obj the callback object to be set.
      */
  void
  setCallbackObject (FaceShiftClient *obj);

  /** \brief Method that checks whether the connection to the FaceShift Studio is active or not. */
  bool
  isConnected ();


signals:
  void
  receivedNetworkUpdates (QSharedPointer<fs::fsMsgTrackingState>);

public slots:
  void
  init ();

  /** \brief Connect to the FaceShift Studio broadcast.
      * \param[in] hostname the address of the host.
      * \param[in] port the host port to connect to.
      */
  void
  connectToServer (QString hostname, int port);


  /** \brief Event that is called when new data is ready to be read. */
  void
  readyRead ();

  /** \brief Event that is called whenever an error occurs in the transmission. */
  void
  error (QAbstractSocket::SocketError);

  /** \brief Event that is called when the state of the socket changes; currently unused, as we use the connected (),
      * error (), and disconnected () events for each specific state.
      */
  void
  stateChanged (QAbstractSocket::SocketState state);

  /** \brief Event that is called when a connection is successfully created. */
  void
  connected ();

  /** \brief Event that is called when the current connection is terminated. */
  void
  disconnected ();

private:
  /** \brief The socket handling the work. */
  QTcpSocket *tcp_socket_;

  /** \brief The size of a message block. */
  quint32 message_block_size_;

  /** \brief Whether there is a new message on the socket, or we are still reading from the same message. */
  bool new_message_;

  /** \brief The object where callbacks are sent for newly received data, after it is de-serialized. */
  FaceShiftClient *callback_object_;

  fsBinaryStream *decoder_;
};

}
